CXX = g++
CXXFLAGS = -std=c++17 -Wall
TARGET = cfd
SRC_DIR = src
BUILD_DIR = build
SRCS = $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SRC_DIR)/Geometry/*.cpp) \
	   $(wildcard $(SRC_DIR)/Array/*.cpp) $(wildcard $(SRC_DIR)/Utilities/*.cpp) \
	   $(wildcard $(SRC_DIR)/Boundary/*.cpp) $(wildcard $(SRC_DIR)/fvm/*.cpp) \
	   $(wildcard $(SRC_DIR)/Matrix/*.cpp) $(wildcard $(SRC_DIR)/StaggeredGrid/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -f $(TARGET)
	rm -rf $(BUILD_DIR)



