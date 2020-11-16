all: gaus
		make cleano
		
gaus: *.o
		g++ *.o -o gaus -lpthread -O3

*.o: *.cpp
		g++ -c *.cpp -lpthread -O3

clean:
		rm -rf *.o gaus
		
cleano:
		rm -rf *.o
