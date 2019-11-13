CC=g++
CFLAGS=-std=c++11
TARGET=Neutronbackground
OBJS=test.C neutron.hxx
LDFLAGS=`root-config --cflags --glibs`

$(TARGET): $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $(OBJS)
