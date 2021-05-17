

def list_to_complex(list):
    real = list[0].replace('(', '').replace(')', '')
    if real[:2] == '--':
        real = real[2:]

    imag = list[2].split('i')[0].replace('(', '')
    if imag[:2] == '--':
        imag = imag[2:]

    opera = float(list[1] + '1')
    return float(real) + opera * float(imag) * 1.j

