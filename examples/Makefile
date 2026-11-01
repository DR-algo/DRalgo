# List of filenames
FILES := \
	01_abelianHiggs.m 03_scalarYukawa.m 05_twoHiggsDoublet.m \
	06_abelianHiggsThermodynamics.m 07_su5GUT.m 10_realScalar.m 11_wessZumino.m \
	12_standardModel.m 13_smCKM.m 14_smZPrime.m 15_smRealSinglet.m \
	16_smTwoSinglets.m 17_smComplexSinglet.m 18_threeHiggsDoublets.m \
	20_higgsTriplet.m 21_leftRightSymmetric.m 23_smSextic.m \
	24_darkSU2HiggsSinglet.m 25_tripletHET.m 26_scotogenic.m

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
