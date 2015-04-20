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


SOURCES  := $(wildcard $(SRCDIR)/*.cpp $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)




$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CFLAGS) $(IPATHS) $< -o $@
	@echo "Compiled "$<" successfully!"


$(TARGET): $(OBJECTS)
	@$(LINKER) $@ $(LFLAGS) $(OBJECTS) 

.PHONY: clean
clean:
	rm $(OBJDIR)/*.o $(TARGET)