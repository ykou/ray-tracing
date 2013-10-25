# sample makefile using Coin3D library

#COIN_DIRECTORY = /usr/local/Coin3d
#COIN_INCLUDE = -I$(COIN_DIRECTORY)/include
#COIN_LIBRARIES = -L$(COIN_DIRECTORY)/lib
#CPP = g++  -R$(COIN_DIRECTORY)/lib
COIN_DIRECTORY = /Library/Frameworks/Inventor.framework/Versions/C
COIN_INCLUDE = -I$(COIN_DIRECTORY)/headers
COIN_LIBRARIES = -L$(COIN_DIRECTORY)/Libraries
CPP = g++  -R$(COIN_DIRECTORY)/Libraries

INCLUDES = $(COIN_INCLUDE) 

# Uncomment to turn on debugging:
CPPFLAGS = -g -DDEBUG

LIBRARIES =  $(COIN_LIBRARIES) -lCoin 

OBJECT_FILES = OSUInventor.o cse681lab4.o
EXEC = rt

.C.o:
	$(CPP) -c $(CPPFLAGS) $(INCLUDES) $<

lab1: $(OBJECT_FILES)
	$(CPP) -o $(EXEC) $(CPPFLAGS) $(OBJECT_FILES) \
        $(LDFLAGS) $(LIBRARIES)

clean:
	rm *.o $(EXEC)
