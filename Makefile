all: gaus
		make cleano
		
gaus: *.o
		g++ *.o -o gaus -lpthread

*.o: *.cpp
		g++ -c *.cpp -lpthread

clean:
		rm -rf *.o gaus
		
cleano:
		rm -rf *.o
