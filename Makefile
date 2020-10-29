
CC=g++

CFLAGS= -pthread -stdlib=libc++ -std=c++11 -m64 -I/Users/ohanabr/Programas/root/include -L/Users/ohanabr/Programas/root/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -lpthread -stdlib=libc++ -lm -ldl -rpath /Users/ohanabr/Programas/root/lib

all:
	$(CC) $(CFLAGS) main.cpp -o main
