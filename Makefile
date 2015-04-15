TARGET = smpor

CXX = g++
LINKER = g++ -o
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

# Modifed Makefile after understanding needed makefile knowledge
SRCDIR = src
OBJDIR = obj


# SOURCES = $(SRCDIR)/*.cpp
# INCLUDES = $(SRCDIR)/*.h
# OBJECTS = $(OBJDIR)/*.o

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
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
	rm *.o $(TARGET)