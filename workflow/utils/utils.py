def get_resource(config, rule, resource) -> int:
	'''
	Function to parse config.yaml to retrieve computational resources for each rule. Returns an int
	'''
	try:
		return config['resources'][rule][resource]
	except KeyError:
		print(f'Failed to get resource for {rule}/{resource}: using default parameters')
		return config['resources']['default'][resource]

