CXX = g++

CXXFLAGS = -mpopcnt -msse4 -mabm -mmmx -mavx -mavx -mavx2 -mfma -O3 
CXXFLAGS2 = -lpthread -mpopcnt -msse4 -mabm -mmmx -mavx -mavx -mavx2 -mfma -O3 

Sufrec: Sufrec_init.o Sufrec_functions_I.o Sufrec_functions_TF.o Sufrec_functions_TC.o Sufrec_main.o
	$(CXX)  Sufrec_init.o Sufrec_functions_I.o Sufrec_functions_TF.o Sufrec_functions_TC.o Sufrec_main.o $(CXXFLAGS2) -o Sufrec 

Sufrec_init.o: Sufrec_init.cpp Sufrec_init.h
	$(CXX) -c Sufrec_init.cpp $(CXXFLAGS) 
Sufrec_functions_I.o: Sufrec_functions_I.cpp Sufrec_init.h Sufrec_functions_I.h
	$(CXX) -c Sufrec_functions_I.cpp $(CXXFLAGS) 
Sufrec_functions_TF.o: Sufrec_functions_TF.cpp Sufrec_init.h Sufrec_functions_TF.h
	$(CXX) -c Sufrec_functions_TF.cpp $(CXXFLAGS2)
Sufrec_functions_TC.o: Sufrec_functions_TC.cpp Sufrec_init.h Sufrec_functions_TC.h
	$(CXX) -c Sufrec_functions_TC.cpp $(CXXFLAGS2)
Sufrec_main.o: Sufrec_main.cpp Sufrec_init.h Sufrec_functions_I.h Sufrec_functions_TF.h Sufrec_functions_TC.h
	$(CXX) -c Sufrec_main.cpp $(CXXFLAGS) 
