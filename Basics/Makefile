all:	debug release

debug:
	g++ settings.cpp main.cpp -g -o test2

release:
	g++ settings.cpp main.cpp -Ofast -DNDEBUG -o test2_release

clean:
	rm test2 test2_release

