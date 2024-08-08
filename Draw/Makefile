# Makefile

CXX := g++
CXXFLAGS := -fPIC -shared `root-config --cflags --libs`
SRC_DIR := src
INC_DIR := interface
LIB_DIR := lib
TARGET := $(LIB_DIR)/libMultiDraw.so
SRC := $(SRC_DIR)/MultiDraw.cc

CWD := $(shell pwd)

INCLUDES := -I$(CWD)/$(INC_DIR)

all: $(LIB_DIR) $(TARGET)

$(TARGET): $(SRC)
	@echo "Compiling $@"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^

$(LIB_DIR):
	@echo "Creating directory $@"
	@mkdir -p $@

clean:
	@echo "Cleaning up..."
	@rm -f $(TARGET)

.PHONY: all clean
