ffmpeg -r 60 -i E/%04d.png -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" movie.mp4

