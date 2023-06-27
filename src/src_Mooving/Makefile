CXX = g++
CXXFLAGS = -mpopcnt -msse4 -mabm -mmmx -mavx -mavx -mavx2 -mfma -O3 

Sufrec: Sufrec_init.o Sufrec_functions_P.o Sufrec_functions_I.o Sufrec_main.o
	$(CXX)  Sufrec_init.o Sufrec_functions_P.o Sufrec_functions_I.o Sufrec_main.o $(CXXFLAGS) -o Sufrec 

Sufrec_init.o: Sufrec_init.cpp Sufrec_init.h
	$(CXX) -c Sufrec_init.cpp $(CXXFLAGS) 
Sufrec_functions_P.o: Sufrec_functions_P.cpp Sufrec_init.h Sufrec_functions_P.h
	$(CXX) -c Sufrec_functions_P.cpp $(CXXFLAGS) 
Sufrec_functions_I.o: Sufrec_functions_I.cpp Sufrec_init.h Sufrec_functions_I.h
	$(CXX) -c Sufrec_functions_I.cpp $(CXXFLAGS) 
Sufrec_main.o: Sufrec_main.cpp Sufrec_init.h Sufrec_functions_P.h Sufrec_functions_I.h 
	$(CXX) -c Sufrec_main.cpp $(CXXFLAGS) 
