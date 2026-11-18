"""Extract the approved result sequence, preserving 60 fps and source time labels."""
from pathlib import Path
import argparse
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ffmpeg', default='ffmpeg')
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists(): p.error('Use a new output filename')
    a.output.parent.mkdir(parents=True, exist_ok=True)
    source = ROOT / 'media/hairdryer_500rpm_landscape_red_blade.mp4'
    cmd = [a.ffmpeg, '-hide_banner', '-loglevel', 'error', '-nostdin', '-threads', '2',
           '-i', str(source), '-an', '-vf', 'trim=start_frame=991,setpts=PTS-STARTPTS,scale=1280:720',
           '-filter_threads', '1', '-c:v', 'libx264', '-preset', 'medium', '-crf', '20',
           '-pix_fmt', 'yuv420p', '-threads', '2', '-movflags', '+faststart', str(a.output)]
    subprocess.run(cmd, check=True)
    print('Summary extracted; source videos and simulation data unchanged.')


if __name__ == '__main__': main()
