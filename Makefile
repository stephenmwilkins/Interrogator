

clean-install:
	pip uninstall flare
	pip uninstall synthobs
	pip install . -r requirements.txt
	rm -rf src

install:
	pip install . -r requirements.txt
	rm -rf src
