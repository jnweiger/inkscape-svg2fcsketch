cd $(dirname $0)
set -x
python ../src/svg2fcsketch.py zigzag.svg > zigzag.fcstd
python ../src/svg2fcsketch.py rect.svg > rect.fcstd
python ../src/svg2fcsketch.py sketch_objects.svg > sketch_objects.fcstd
python ../src/svg2fcsketch.py arcs.svg > arcs.fcstd
python ../src/svg2fcsketch.py rects.svg -i rect1384_origin > rect1384.fcstd
python ../src/svg2fcsketch.py paths_study2.svg > paths_study2.fcstd
