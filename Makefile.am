# GNU Automake config
SUBDIRS = data gctpc/source src test cvmh_test
INCLUDES = $(default_includes)

.PHONY = run_api run_unit run_accept

run_api:
	cd tests;$(MAKE) run_api

run_unit:
	cd test;$(MAKE) run_unit

run_accept:
	cd test;$(MAKE) run_accept
