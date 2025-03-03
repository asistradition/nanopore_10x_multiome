RC_SIMPLE = {
    ord('A'): 'T',
    ord('T'): 'A',
    ord('G'): 'C',
    ord('C'): 'G',
    ord('N'): 'N',
    ord('a'): 't',
    ord('t'): 'a',
    ord('g'): 'c',
    ord('c'): 'g',
    ord('n'): 'n'
}


def RC(x):
    return x.translate(RC_SIMPLE)[::-1]


def REV(x):
    return x[::-1]
