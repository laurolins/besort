nothing:

waf:
	wget http://ftp.waf.io/pub/release/waf-1.7.16
	mv waf-1.7.16 waf
	python waf --version
	chmod +x waf

update:
	./waf build
	./waf install

all:
	./waf distclean
	./waf configure --prefix=${HOME}/local-dev
	./waf build
	./waf install

verbose:
	./waf distclean
	./waf configure --prefix=${HOME}/local-dev
	./waf -vvv build
	./waf install

test:
	TILETREE_ROOT="${HOME}/projects/tiletree" ./build/src/tiletree
