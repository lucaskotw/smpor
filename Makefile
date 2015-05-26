TARGET = smpor

CC = gcc
CXX = g++
LINKER = g++ -o
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall -framework OpenGL -lglfw3 $(DEBUG)
IPATHS = -I/usr/local/Cellar/eigen/3.2.4/include/eigen3/ -I/usr/local/include -I/opt/X11/include -I/usr/local/Cellar/devil/1.7.8_1/include/
LPATHS = -L/usr/local/lib -L/opt/X11/lib

# Modifed Makefile after understanding needed makefile knowledge
SRCDIR = src
OBJDIR = obj


CXXSOURCES  := $(wildcard $(SRCDIR)/*.cpp)
CSOURCES  := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(SRCDIR)/*.h)
CXXOBJECTS  := $(CXXSOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
COBJECTS  := $(CSOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)



$(CXXOBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CFLAGS) $(IPATHS) $< -o $@
	@echo "Compiled "$<" successfully!"


$(COBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) $(IPATHS) $< -o $@
	@echo "Compiled "$<" successfully!"



$(TARGET): $(CXXOBJECTS) $(COBJECTS)
	@$(LINKER) $@ $(LFLAGS) $(CXXOBJECTS) $(COBJECTS)

	

.PHONY: clean
clean:
	rm $(OBJDIR)/*.o $(TARGET)