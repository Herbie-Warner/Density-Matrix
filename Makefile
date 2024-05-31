CXX := g++
CXXFLAGS := -Iinclude -IC:\Users\herbi\eigen-3.4.0
SRCDIR := src
OBJDIR := obj
BINDIR := bin

METHODS_SRCS := $(wildcard $(SRCDIR)/methods/*.cpp)
MAIN_SRC := $(wildcard $(SRCDIR)/main/*.cpp)

SRCS := $(METHODS_SRCS) $(MAIN_SRC)
OBJS := $(SRCS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
OBJS := $(OBJS:/=\)
TARGET := optimisation.exe

all: build $(TARGET)

build:
	@if not exist "$(OBJDIR)\methods" mkdir "$(OBJDIR)\methods"
	@if not exist "$(OBJDIR)\main" mkdir "$(OBJDIR)\main"
	@if not exist "$(BINDIR)" mkdir "$(BINDIR)"

$(TARGET): $(OBJS)
	$(CXX) $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if not exist "$(dir $@)" mkdir "$(dir $@)"
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	@if exist "$(OBJDIR)" rmdir /s /q "$(OBJDIR)"
	@if exist "$(BINDIR)" rmdir /s /q "$(BINDIR)"

.PHONY: all build clean
