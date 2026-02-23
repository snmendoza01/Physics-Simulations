import imageio
import os

frames_dir = "frames"

frame_files = sorted([
    f for f in os.listdir(frames_dir)
    if f.endswith(".png")
])

with imageio.get_writer("animation.mp4", fps=15) as writer:
    for filename in frame_files:
        filepath = os.path.join(frames_dir, filename)
        image = imageio.imread(filepath)
        writer.append_data(image)