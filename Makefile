Neutronbackgrond: test.C neutron.hxx
	g++ -o Neutronbackground -std=c++11 `root-config --cflags --glibs` test.C
