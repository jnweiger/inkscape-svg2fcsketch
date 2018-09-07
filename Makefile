#

EXTNAME=svg2fcsketch
DISTNAME=inkscape-$(EXTNAME)
VERS=$$(python src/svg2fcsketch.py --version | sed -n -e 's@.*\.py\s*@@p')

all: clean build check

build: $(EXTNAME).py $(EXTNAME).inx

$(EXTNAME).inx:
	cp src/$(EXTNAME).inx .
	@grep v$(VERS) $(EXTNAME).inx || echo "WARNING: Version number v$(VERS) not seen in $(EXTNAME).inx"

$(EXTNAME).py:
	@sed >  $@ -e '/INLINE_BLOCK_START/,$$d' < src/$(EXTNAME).py
	sed >> $@ -e '/if __name__ ==/,$$d' < inksvg/src/inksvg.py
	@sed >> $@ -e '1,/INLINE_BLOCK_END/d' < src/$(EXTNAME).py

# selftest
check:
	test/RUN.sh

clean:
	rm -f $(EXTNAME).py *.inx
	rm -f *.orig */*.orig
	rm -rf distribute/$(DISTNAME)
	rm -rf distribute/deb/files

