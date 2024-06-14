# List of filenames
FILES := \
	2hdm.m 2xsm.m 3hdm.m ah-thermo.m ah.m cxsm.m htm.m LRSymmetric.m \
	RealScalar.m ScalYukawa.m sm.m \
	smCKM.m smDim6.m smZp.m SU5.m \
	tripletHET.m WessZumino.m xsm.m \
	dark-su2-higgs-singlet.m

# Target: test
test: make

# Target: make
make: create_log_folder $(foreach file,$(FILES),run_wolframscript_$(file))

# Target: create_log_folder
create_log_folder:
	@mkdir -p log

# Rule for each filename
define run_wolframscript_template
.PHONY: run_wolframscript_$(1)
run_wolframscript_$(1): create_log_folder
	@echo "Running wolframscript on $(1)"
	@wolframscript -file $(1) >> log/$(1).log
endef

# Apply the template rule for each file
$(foreach file,$(FILES),$(eval $(call run_wolframscript_template,$(file))))

clean:
	rm -rf log/*

.PHONY: test clean
