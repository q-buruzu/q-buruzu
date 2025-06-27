CXX = g++
CXXFLAGS = -Wall
BUILD_NAME = qburuzu
BUILD_DIR = build
SRC_DIR = src

SRCS := $(wildcard $(SRC_DIR)/*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))

all: $(BUILD_NAME)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(BUILD_DIR)/$@

run:
	./$(BUILD_DIR)/$(BUILD_NAME)

clean:
	rm -rf $(BUILD_DIR)
	mkdir $(BUILD_DIR)
