MY_NOSE_FLAGS?=-v -s

wheel:
	# This is basically a syntax check.
	python setup.py bdist_wheel
doctest:
	py.test ${MY_TEST_FLAGS} --doctest-modules svkits/
pylint:
	pylint --errors-only svkits/
unit-test:
	# These use the local tree. No installation required. Extremely fast.
	nosetests ${MY_NOSE_FLAGS} utests.py
cram:
	cd tests/cram && cram tests/cram/*.t && cd ..
autofmt:
	find svkits -type f -name '*.py' | xargs autoflake --in-place --remove-unused-variables --expand-star-imports
	find svkits -type f -name '*.py' | xargs autopep8  --in-place --max-line-length 120
vulture:
	vulture svkits/
clean:
	find . -type f -name '*.pyc' |xargs rm -f
