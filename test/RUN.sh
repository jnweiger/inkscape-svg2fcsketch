cd $(dirname $0)
set -x
python ../src/svg2fcsketch.py zigzag.svg
python ../src/svg2fcsketch.py rect.svg
python ../src/svg2fcsketch.py sketch_objects.svg
python ../src/svg2fcsketch.py arcs.svg
python ../src/svg2fcsketch.py rects.svg -i rect1384_origin
