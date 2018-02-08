MY_NOSE_FLAGS?=-v -s

wheel:
	# This is basically a syntax check.
	python setup.py bdist_wheel
pylint:
	pylint --errors-only svkits/
unit-test:
	# These use the local tree. No installation required. Extremely fast.
	nosetests ${MY_NOSE_FLAGS} utests.py
cram:
	cd tests/cram && cram tests/cram/*.t && cd ..
