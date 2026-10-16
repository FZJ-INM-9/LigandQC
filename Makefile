update_sphinx_documentation:
	rm -rf docs/source/api && mkdir -p docs/source/api
	sphinx-apidoc -o docs/source/api ligandqc --separate --force --module-first --maxdepth 1
	rm -rf docs/build && cd docs && make html

open_sphinx_documentation:
	open docs/build/html/index.html

live_sphinx_documentation:
	cd docs && make livehtml