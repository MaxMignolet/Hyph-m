!#/bin/sh

cd dir_viz

ffmpeg -r 10 -i image_20_%4d_top.png  -c:v libx264 -vf "scale=-2:'min(1080,ih)',fps=10" -pix_fmt yuv420p anim_20_top.mp4
ffmpeg -r 10 -i image_20_%4d_front.png  -c:v libx264 -vf "scale=-2:'min(1080,ih)',fps=10" -pix_fmt yuv420p anim_20_front.mp4

ffmpeg -r 10 -i image_21_%4d_top.png  -c:v libx264 -vf "scale=-2:'min(1080,ih)',fps=10" -pix_fmt yuv420p anim_21_top.mp4
ffmpeg -r 10 -i image_21_%4d_front.png  -c:v libx264 -vf "scale=-2:'min(1080,ih)',fps=10" -pix_fmt yuv420p anim_21_front.mp4

ffmpeg -r 10 -i image_22_%4d_top.png  -c:v libx264 -vf "scale=-2:'min(1080,ih)',fps=10" -pix_fmt yuv420p anim_22_top.mp4
ffmpeg -r 10 -i image_22_%4d_front.png  -c:v libx264 -vf "scale=-2:'min(1080,ih)',fps=10" -pix_fmt yuv420p anim_22_front.mp4

ffmpeg -r 10 -i image_23_%4d_top.png  -c:v libx264 -vf "scale=-2:'min(1080,ih)',fps=10" -pix_fmt yuv420p anim_23_top.mp4
ffmpeg -r 10 -i image_23_%4d_front.png  -c:v libx264 -vf "scale=-2:'min(1080,ih)',fps=10" -pix_fmt yuv420p anim_23_front.mp4

cd ..
