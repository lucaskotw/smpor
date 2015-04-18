TARGET = smpor

CC = gcc
CXX = g++
LINKER = g++ -o
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall -framework OpenGL -lglfw3 $(DEBUG)
IPATHS = -I/usr/local/Cellar/eigen/3.2.4/include/eigen3/ -I/usr/local/include -I/opt/X11/include
LPATHS = -L/usr/local/lib -L/opt/X11/lib

# Modifed Makefile after understanding needed makefile knowledge
SRCDIR = src
OBJDIR = obj


# SOURCES = $(SRCDIR)/*.cpp
# INCLUDES = $(SRCDIR)/*.h
# OBJECTS = $(OBJDIR)/*.o

SOURCES  := $(wildcard $(SRCDIR)/*.cpp $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

# $(TARGET): $(OBJS)



$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	@$(CXX) $(CFLAGS) $< -o $@
	@echo "Compiled "$<" successfully!"
# graph.o: src/graph.h src/graph.cpp 
# 	$(CXX) $(CFLAGS) src/graph.cpp

# main.o: src/graph.h src/main.cpp
# 	$(CXX) $(CFLAGS) src/main.cpp

$(TARGET): $(OBJECTS)
	@$(LINKER) $@ $(LFLAGS) $(OBJECTS) 
	# $(CXX) graph.o main.o $(LFLAGS) -o smpor

.PHONY: clean
clean:
	rm $(OBJDIR)/*.o $(TARGET)