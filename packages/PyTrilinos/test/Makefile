# Output redirect macro
TEE := | tee -a

# Target for testing without doing a diff
test_no_diff:
	$(MAKE) TEST_NO_DIFF=1 test

# Diff testing macro
TEST_DIFF := ./numSavvyDiff $(TEST_OUT) $(TEST_OUT).gold

# Target for testing with doing a diff
test_diff:
# Put if test here so can do "make test_no_diff" and no diffing is done
# This is a good way to check if all tests run to completion without a bad
# exit status
ifeq ($(TEST_NO_DIFF),)
	@echo
	@echo --------------------
	@echo $(TEST_DIFF)
	@     $(TEST_DIFF);            \
	if [ "$$?" = "0" ]; then       \
	  echo "No differences found"; \
	  echo --------------------;   \
	  echo;                        \
	else                           \
	  echo "Differences found";    \
	  echo --------------------;   \
	  echo;                        \
	  exit 1;                      \
	fi
endif

test: testNonContiguous testEpetra testSerialDense testParameter testStatusTest \
      testExample

testNonContiguous:
	rm -fr $@.out
	touch $@.out
	@echo "----" $(TEE) $@.out
	@echo "Test non-contiguous to contiguous NumPy array translation" $(TEE) $@.out
	@echo "----" $(TEE) $@.out
	@./$@.py 2>&1 $(TEE) $@.out
	$(MAKE) TEST_OUT=$@.out test_diff

testEpetra:
	rm -fr $@.out
	touch $@.out
	@echo "----" $(TEE) $@.out
	@echo "Test Epetra interface" $(TEE) $@.out
	@echo "----" $(TEE) $@.out
	@./$@.py 2>&1 $(TEE) $@.out
	$(MAKE) TEST_OUT=$@.out test_diff

testSerialDense:
	rm -fr $@.out
	touch $@.out
	@echo "----" $(TEE) $@.out
	@echo "Test the SerialDense interface" $(TEE) $@.out
	@echo "----" $(TEE) $@.out
	@./$@.py 2>&1 $(TEE) $@.out
	$(MAKE) TEST_OUT=$@.out test_diff

testParameter:
	rm -fr $@.out
	touch $@.out
	@echo "----" $(TEE) $@.out
	@echo "Test Parameter.List interface" $(TEE) $@.out
	@echo "----" $(TEE) $@.out
	@./$@.py 2>&1 $(TEE) $@.out
	$(MAKE) TEST_OUT=$@.out test_diff

testStatusTest:
	rm -fr $@.out
	touch $@.out
	@echo "----" $(TEE) $@.out
	@echo "Test the StatusTest interface" $(TEE) $@.out
	@echo "----" $(TEE) $@.out
	@./$@.py 2>&1 $(TEE) $@.out
	$(MAKE) TEST_OUT=$@.out test_diff

testExample:
	rm -fr $@.out
	touch $@.out
	@echo "----" $(TEE) $@.out
	@echo "Test an example NOX_Epetra nonlinear problem" $(TEE) $@.out
	@echo "----" $(TEE) $@.out
	@./$@.py 2>&1 $(TEE) $@.out
	@grep -v time $@.out > TEMP
	@mv TEMP $@.out
	$(MAKE) TEST_OUT=$@.out test_diff
