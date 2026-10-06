.PHONY: doctor test notebook probe
doctor:
	./scripts/setup/doctor.sh
test:
	python -m pytest -q
notebook:
	jupyter lab
probe:
	./scripts/setup/probe_server.sh
