ffmpeg -r 1 -f image2 -s 600x600 -i ./pict/img_%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p ./movie.mp4
