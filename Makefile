.PHONY: clean clean-test clean-pyc clean-build docs help
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python -c "$$BROWSER_PYSCRIPT"

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr .tox/
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache

lint: ## check style with flake8
	python3 -m flake8 mps tests

type: ## Run mypy
	python3 -m mypy mps tests

black: ## Run mypy
	python3 -m black mps tests

test: ## run tests quickly with the default Python
	python3 -m pytest --cov=mps tests

docs: ## generate Sphinx HTML documentation, including API docs
	rm -f docs/mps.rst
	rm -f docs/modules.rst
	sphinx-apidoc -o docs/ mps
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

servedocs: docs ## compile the docs watching for changes
	watchmedo shell-command -p '*.rst' -c '$(MAKE) -C docs html' -R -D .

dist: clean ## builds source and wheel package
	python setup.py sdist
	python setup.py bdist_wheel
	ls -l dist

install: clean ## install the package to the active Python's site-packages
	python -m pip install .

dev: clean ## Just need to make sure that libfiles remains
	python -m pip install --upgrade pip
	python3 -m pip install -r requirements_dev.txt
	python3 -m pip install -e ".[all]"
	pre-commit install

installer: clean
	python -m pip install "matplotlib==3.0.3"
	python -m pip install -e ".[all]"
	python -m pip install pyinstaller
	pyinstaller -F mps/__main__.py -n mps

installer-windows: clean
	python -m pip install pipwin
	pipwin install "matplotlib==3.0.3"
	pipwin install -e ".[all]"
	pipwin install pyinstaller
	pyinstaller -F mps/__main__.py -n mps
