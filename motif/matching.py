def find_motif_positions_dfa(seq, dfa, motif_len, alphabet):
    positions = []
    state = dfa.start
    for j, ch in enumerate(seq.upper()):
        if ch not in alphabet:
            state = dfa.start
            continue
        key = (state, ch)
        if key in dfa.delta:
            state = dfa.delta[key]
        else:
            state = dfa.start
            key2 = (state, ch)
            if key2 in dfa.delta:
                state = dfa.delta[key2]
        if state in dfa.accepts:
            positions.append(j - motif_len + 1)
    return [p for p in positions if p >= 0]
