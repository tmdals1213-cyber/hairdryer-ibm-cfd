"""LICA staggered-field reader and nonuniform-grid Q/lambda2 diagnostics.

Coordinates and velocities remain in solver units. For Uref=Lref=1, Q and
lambda2 are nondimensional; multiply by (Uref/Lref)^2 for SI units.
"""
from pathlib import Path
import struct
import numpy as np


def read_grid(path):
    tok=Path(path).read_text().split();dims=tuple(map(int,tok[:3]));pos=6;axes=[]
    for n in dims:
        a=np.array(tok[pos:pos+n],float);pos+=n
        if len(a)!=n or not np.all(np.diff(a)>0):raise ValueError('Invalid grid coordinates')
        axes.append(a)
    return axes


def read_records(path):
    records=[]
    with Path(path).open('rb') as f:
        while marker:=f.read(4):
            if len(marker)!=4:raise ValueError('Truncated record marker')
            n=struct.unpack('<i',marker)[0]
            if n<0 or n>4_000_000_000:raise ValueError('Invalid record length')
            data=f.read(n)
            if len(data)!=n or f.read(4)!=marker:raise ValueError('Corrupt Fortran record')
            records.append(data)
    return records


def read_field(path,axes,include_gradient=False):
    x,y,z=axes;nx,ny,nz=map(len,axes);a=read_records(path)
    if len(a)!=8:raise ValueError('Expected stationary LICA eight-record field')
    n1,n2,n3,re=struct.unpack('<qqqd',a[0])
    if (n1,n2,n3)!=(nx,ny,nz):raise ValueError('Grid/field dimensions mismatch')
    step,segment_step,time,dt=struct.unpack('<qqdd',a[1])
    u=np.frombuffer(a[4],'<f8').reshape((nx,ny+1,nz+1),order='F')
    v=np.frombuffer(a[5],'<f8').reshape((nx+1,ny,nz+1),order='F')
    w=np.frombuffer(a[6],'<f8').reshape((nx+1,ny+1,nz),order='F')
    velocity=np.stack([(u[:-1,1:ny,1:nz]+u[1:,1:ny,1:nz])/2,
        (v[1:nx,:-1,1:nz]+v[1:nx,1:,1:nz])/2,
        (w[1:nx,1:ny,:-1]+w[1:nx,1:ny,1:])/2],axis=-1)
    pressure=np.frombuffer(a[7],'<f8').reshape((nx-1,ny-1,nz-1),order='F')
    divergence=np.diff(u[:,1:ny,1:nz],axis=0)/np.diff(x)[:,None,None]+np.diff(v[1:nx,:,1:nz],axis=1)/np.diff(y)[None,:,None]+np.diff(w[1:nx,1:ny,:],axis=2)/np.diff(z)[None,None,:]
    if not np.isfinite(velocity).all() or not np.isfinite(pressure).all():raise ValueError('Nonfinite flow field')
    inlet=float(np.sum(u[0,1:ny,1:nz]*np.diff(y)[:,None]*np.diff(z)[None,:]))
    net=float(np.sum((u[-1,1:ny,1:nz]-u[0,1:ny,1:nz])*np.diff(y)[:,None]*np.diff(z)[None,:])+np.sum((v[1:nx,-1,1:nz]-v[1:nx,0,1:nz])*np.diff(x)[:,None]*np.diff(z)[None,:])+np.sum((w[1:nx,1:ny,-1]-w[1:nx,1:ny,0])*np.diff(x)[:,None]*np.diff(y)[None,:]))
    result={'velocity':velocity,'pressure':pressure,'finite_volume_divergence':divergence,'Re':re,'time':time,'dt':dt,'step':step,'segment_step':segment_step,'inlet_flux':inlet,'net_outward_boundary_flux':net,'relative_boundary_flux_imbalance':abs(net)/max(abs(inlet),1e-30),'centres':[(a[1:]+a[:-1])/2 for a in axes]}
    if include_gradient:
        result['velocity_gradient']=staggered_gradient(velocity,result['centres'],[u[:,1:ny,1:nz],v[1:nx,:,1:nz],w[1:nx,1:ny,:]],axes)
    return result


def staggered_gradient(velocity,centres,face_velocity,face_coordinates):
    """Cell-centred gradient retaining the native MAC normal face differences.

    Off-diagonal terms use second-order derivatives of face-averaged velocity.
    Diagonal terms are direct opposing-face differences / cell width, so the
    tensor trace is exactly the solver's finite-volume divergence. On stretched
    meshes differentiating the already averaged velocity need not preserve it.
    """
    G=np.stack([np.stack(np.gradient(velocity[...,i],*centres,edge_order=2),axis=-1) for i in range(3)],axis=-2)
    for i,(values,nodes) in enumerate(zip(face_velocity,face_coordinates)):
        shape=[1,1,1];shape[i]=len(nodes)-1
        G[...,i,i]=np.diff(values,axis=i)/np.diff(nodes).reshape(shape)
    return G


def vortex_criteria(velocity,centres,gradient=None):
    """G_ij=du_i/dx_j, S=(G+G.T)/2, Omega=(G-G.T)/2.

    Q=0.5*(||Omega||_F^2-||S||_F^2).
    lambda2 is the MIDDLE eigenvalue of symmetric S@S+Omega@Omega.
    np.gradient uses physical nonuniform coordinates, second-order formulas.
    """
    v=np.asarray(velocity,float)
    if v.shape!=tuple(map(len,centres))+(3,):raise ValueError('Velocity/coordinate mismatch')
    G=np.stack([np.stack(np.gradient(v[...,i],*centres,edge_order=2),axis=-1) for i in range(3)],axis=-2) if gradient is None else np.asarray(gradient,float)
    if G.shape!=v.shape[:-1]+(3,3):raise ValueError('Gradient tensor shape mismatch')
    S=.5*(G+G.swapaxes(-1,-2));Omega=.5*(G-G.swapaxes(-1,-2))
    Q=.5*(np.sum(Omega*Omega,axis=(-2,-1))-np.sum(S*S,axis=(-2,-1)))
    A=S@S+Omega@Omega;lambda2=np.linalg.eigvalsh(.5*(A+A.swapaxes(-1,-2)))[...,1]
    curl=np.stack([G[...,2,1]-G[...,1,2],G[...,0,2]-G[...,2,0],G[...,1,0]-G[...,0,1]],axis=-1)
    return {'Q':Q,'lambda2':lambda2,'vorticity':curl,'gradient_divergence':np.trace(G,axis1=-2,axis2=-1)}


def stencil_valid(signed_distance,wall_guard=0.):
    """Only all-fluid central stencils; discard outer derivative boundary nodes."""
    fluid=np.asarray(signed_distance)>0;valid=fluid&(signed_distance>wall_guard)
    for axis in range(3):
        valid &= np.roll(fluid,1,axis)&np.roll(fluid,-1,axis)
        for index in [0,-1]:
            sl=[slice(None)]*3;sl[axis]=index;valid[tuple(sl)]=False
    return valid
