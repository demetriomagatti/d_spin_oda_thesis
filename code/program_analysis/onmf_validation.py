def validate_onmf(model):
    print(type(model.onmf_decomposition))
    print(model.onmf_decomposition.get_params())
    print("components_:", model.onmf_decomposition.components_.shape)