html:
	sphinx-build -b html ./source/ ./generated_docs/docs -W --keep-going
	python3 ./source/edit_button_handler.py

clean:
	rm -rf ./generated_docs/
