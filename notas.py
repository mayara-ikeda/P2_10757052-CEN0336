def calcular_media():
    try:
        # input das notas entre 0 e 10
        notas_input = input("Digite todas as notas separadas por espaço (valores entre 0 e 10, use ponto ou vírgula como separador decimal): ")
        # Converte vírgulas para pontos e cria uma lista de floats
        notas = [float(nota.replace(',', '.')) for nota in notas_input.split()]
        
        if all(0 <= nota <= 10 for nota in notas):  # Valida que todas as notas estão no intervalo entre 0 e 10
            total = sum(notas)  # Calcula o total das notas
            quantidade = len(notas)  # Quantidade total de notas para calcular a média (permite input de quantas notas for necessário, seja mais ou menos que 10 notas)
            media = total / quantidade  # Calcula a média
            print(f"A média da disciplina é: {media:.2f}")  # Exibe a média formatada
        else:
            print("Erro: Todas as notas devem estar no intervalo de 0 a 10. Tente novamente.") # a condição cai neste else caso alguma nota esteja fora do range de 0 a 10
    except ValueError:
        print("Erro: Certifique-se de inserir apenas números separados por espaço.") # A condição cai neste Except quando os números estão separados de forma incorreta (ponto e vírgula ou barra por ex)

# Executa a função principal, evitando sair da main
if __name__ == "__main__":
    calcular_media()
