LIB_DIR = lib

default: icsi2022

icsi2022: setup.py icsi2022.pyx $(LIB_DIR)/icsi2022.a
	python setup.py build_ext --inplace && rm -f icsi2022.c && rm -Rf build

$(LIB_DIR)/icsi2022.a:
	make -C $(LIB_DIR) libicsi2022.a

clean:
	rm *.so