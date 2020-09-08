Shift = 16
Beta = 4135
Alpha = Base - 3 * Base

def guess_by_b(N, seq_a, seq_std):
    
    begin_num = 0
    # first 32 guess is always 0
    for i in range(0, 32):
        correct_count += 1 == seq_std[i]
    for i in range(1, N):


def info_by_a(N, seq_ans):

    the_slice_num = 2 ** 20
    last_trunk = 2 ** 10
    assert N = the_slice_num

    seq_a = [False] * N
    seq_b = [False] * N
    trunk_indicator = 0
    trunk_base_shift = last_trunk
    for i in range(N - last_trunk, N):
        seq_a[i] = seq_ans[i]
        seq_b[i] = seq_ans[i]
        trunk_indicator = trunk_indicator * 2 + seq_ans[i]
    N -= last_trunk 
    for i in range(0, 32):
        index = N - 32 + i
        seq_a[index] = (last_trunk >> i) & 1
    
    N -= 32
    cases = []
    
    while True:
        multi = trunk_indicator * 2 ** (trunk_base_shift)
    
    

