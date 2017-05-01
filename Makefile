#!/usr/bin/make -f

workflow_diagram.png: workflow_diagram.dot
	dot -Tpng $< -o $@
