OBJS = Main.o General.o Trees_Edges.o Manipulate_Matrix.o Matrix_Scan.o Parse_Data.o User_Interface.o Parse_Nexus.o User_Tree.o GAoptimize.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c -m64 -O3 -funroll-loops $(DEBUG)
LFLAGS = -Wall -m64 $(DEBUG)

Decisivator: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o Decisivator_serial
	
Main.o: Main.cpp
	$(CC) $(CFLAGS) Main.cpp
	
General.o: General.cpp General.h
	$(CC) $(CFLAGS) General.cpp
	
Trees_Edges.o: Trees_Edges.cpp Trees_Edges.h
	$(CC) $(CFLAGS) Trees_Edges.cpp
	
Manipulate_Matrix.o: Manipulate_Matrix.cpp Manipulate_Matrix.h
	$(CC) $(CFLAGS) Manipulate_Matrix.cpp
	
Matrix_Scan.o: Matrix_Scan.cpp Matrix_Scan.h
	$(CC) $(CFLAGS) Matrix_Scan.cpp
	
Parse_Data.o: Parse_Data.cpp Parse_Data.h
	$(CC) $(CFLAGS) Parse_Data.cpp
	
User_Interface.o: User_Interface.cpp User_Interface.h
	$(CC) $(CFLAGS) User_Interface.cpp
	
Parse_Nexus.o: Parse_Nexus.cpp Parse_Nexus.h
	$(CC) $(CFLAGS) Parse_Nexus.cpp

User_Tree.o: User_Tree.cpp User_Tree.h
	$(CC) $(CFLAGS) User_Tree.cpp

GAoptimize.o: GAoptimize.cpp GAoptimize.h
	$(CC) $(CFLAGS) GAoptimize.cpp

clean:
	rm -rf *.o Decisivator_serial