# List of filenames
FILES := \
	01_abelian_higgs.m 03_scalar_yukawa.m 05_two_higgs_doublet.m \
	06_ah_thermodynamics.m 07_su5_gut.m 10_real_scalar.m 11_wess_zumino.m \
	12_standard_model.m 13_sm_ckm.m 14_sm_z_prime.m 15_sm_real_singlet.m \
	16_sm_two_singlets.m 17_sm_complex_singlet.m 18_three_higgs_doublets.m \
	20_higgs_triplet.m 21_left_right_symmetric.m 23_sm_sextic.m \
	24_dark_su2_higgs_singlet.m 25_triplet_het.m 26_scotogenic.m

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
