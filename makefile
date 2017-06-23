MY_NOSE_FLAGS?=-v -s

wheel:
	# This is basically a syntax check.
	python setup.py bdist_wheel

unit-test:
	# These use the local tree. No installation required. Extremely fast.
	nosetests ${MY_NOSE_FLAGS} utests.py

pylint:
	pylint --errors-only svkits/*.py
