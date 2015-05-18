% Generate some data
N = 100; % Number of sequences
L = 25;  % Sequence length
[seq, labels] = generate_data(N, L);
labels_est = sphere_clustering(seq, 2, 1, 'linear');
