# variáveis repetidas
CC=gcc
CCFLAGS=-Wall -Wno-unused-result
CCLIBS=-lm
BINFILE=ep1

# compila e executa, regra padrão
executar:
	make limpar
	make compilar && ./$(BINFILE)

# compila e debuga com gdb
debugar_gdb:
	make limpar
	make compilar_debug && gdb ./$(BINFILE)

# compila e debuga com valgrind (memcheck)
debugar_valgrind:
	make limpar
	make compilar_debug && valgrind ./$(BINFILE)

# remove os binários
limpar:
	rm -f $(BINFILE)

# compila para execução normal
compilar:
	$(CC) *.c $(CCFLAGS) -O3 -o $(BINFILE) $(CCLIBS)

# compila para debug
compilar_debug:
	$(CC) *.c $(CCFLAGS) -g -o $(BINFILE) $(CCLIBS)
