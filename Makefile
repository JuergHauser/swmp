NOTEBOOKS := $(wildcard examples/*.py)
RENDERED  := $(patsubst examples/%.py,examples/rendered/%.html,$(NOTEBOOKS))

.PHONY: render-notebooks clean-notebooks

render-notebooks: $(RENDERED)

examples/rendered/%.html: examples/%.py | examples/rendered
	marimo export html $< -o $@ --include-code --force

examples/rendered:
	mkdir -p $@

clean-notebooks:
	rm -f $(RENDERED)
