CXXFLAGS = -Wall -std=c++14 -I. -O2
BASE_SRCS = $(wildcard sampling/*.cpp)

all: generate_samples shuffle_indices

release: CXXFLAGS += -O3
release: generate_samples shuffle_indices

generate_samples: generate_samples.cpp $(BASE_SRCS)
	g++ $(CXXFLAGS) -o generate_samples generate_samples.cpp $(BASE_SRCS)

shuffle_indices: shuffle_indices.cpp $(BASE_SRCS)
	g++ $(CXXFLAGS) -o shuffle_indices shuffle_indices.cpp $(BASE_SRCS)

clean:
	rm -Rf generate_samples shuffle_indices