.PHONY: all nfa_engine generator copy

all: nfa_engine generator copy

nfa_engine:
	cd nfa_engine && $(MAKE)
	
generator:
	cd generator && $(MAKE)

copy:
	cp nfa_engine/nfa_engine bin/	
	cp generator/regex_memory bin/

clean:
	rm -f bin/nfa_engine bin/regex_memory
	cd generator && $(MAKE) clean
	cd nfa_engine && $(MAKE) clean
