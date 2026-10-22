
TUTORIALS := 1 2 3 4 5 6


update_sphinx_documentation:
	rm -rf docs/source/api && mkdir -p docs/source/api
	sphinx-apidoc -o docs/source/api ligandqc --separate --force --module-first --maxdepth 1
	rm -rf docs/build && cd docs && make html

open_sphinx_documentation:
	open docs/build/html/index.html

live_sphinx_documentation:
	cd docs && make livehtml

run_tests:
	python -m pytest tests --maxfail=1 --disable-warnings

lint:
	ruff check ligandqc --config pyproject.toml

lint_fix:
	ruff check ligandqc --config pyproject.toml --fix

run_tests_complete:
	python -m pytest tests --disable-warnings

format:
	black ligandqc --config pyproject.toml

all_tests: format lint run_tests_complete

run_tutorial_%:
	cd docs/source/_static/tutorial_files/tut$* && bash run.sh

run_tutorials: $(addprefix run_tutorial_, $(TUTORIALS))

run_everything: all_tests run_tutorials update_sphinx_documentation