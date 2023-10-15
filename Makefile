build:
	python3 setup.py build_ext --inplace

run:
	python3 main.py

build_and_run: build run