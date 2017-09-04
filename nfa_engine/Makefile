.PHONY: release real

CUDA_OBJ = cuda_wrappers gpu_functions host_functions nfa_engine burst
HOST_OBJ = cuda_allocator config_options globals half_trie sparse transition_graph state_vector utils
COMMON_HEADERS = common.h

NVCC=nvcc
#SM=sm_35
#SM=sm_52
SM=sm_61
CXX=g++
CC=gcc

######CUDA_INCLUDE=-I/usr/local/cuda/include
CUDA_INCLUDE=-I/usr/local/cuda-8.0/include

CXXFLAGS+=$(CUDA_INCLUDE) -Wno-deprecated
NVCCFLAGS+=$(CUDA_INCLUDE) -Xptxas -v -arch ${SM} --compiler-options -Wno-deprecated -lineinfo -DTEXTURE_MEM_USE

export DYN_LIB=-lpcap -lboost_regex -lpcre

NVCCFLAGS_REL = $(NVCCFLAGS) -O4

CXXFLAGS_REL = $(CXXFLAGS) -O4

release:
	$(MAKE) -e real NVCCFLAGS="$(NVCCFLAGS_REL)" CXXFLAGS="$(CXXFLAGS_REL)"

real: nfa_engine

$(addsuffix .o, $(HOST_OBJ)) $(addsuffix .o, $(CUDA_OBJ)) : $(COMMON_HEADERS)

$(addsuffix .o, $(HOST_OBJ)) : $(addsuffix .cpp, $(basename $@)) $(addsuffix .h, $(basename $@))

$(addsuffix .o, $(CUDA_OBJ)) : $(addsuffix .cu, $(basename $@)) $(addsuffix .h, $(basename $@))
	${NVCC} $(NVCCFLAGS) -c -o $(addsuffix .o, $(basename $@)) $(addsuffix .cu, $(basename $@))
	
nfa_engine: $(addsuffix .o, $(HOST_OBJ)) $(addsuffix .o, $(CUDA_OBJ))
	${NVCC} $(NVCCFLAGS) -o nfa_engine $(addsuffix .o, $(HOST_OBJ)) $(addsuffix .o, $(CUDA_OBJ)) ${DYN_LIB}

clean:
	rm -f *.o nfa_engine
