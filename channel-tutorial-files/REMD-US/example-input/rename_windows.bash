for i in {0..18}; do
	for ext in vel xsc coor; do
		ln -s window${i}_min.${ext} AmtB-${i}.restart.${ext}
	done
done
